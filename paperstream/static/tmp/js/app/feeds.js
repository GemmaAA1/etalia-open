define(
    ['jquery', 'app/ui/detail', 'app/ui/layout', 'app/util/utils', 'app/templates', 'endless', 'bootstrap'],
    function($, Detail, Layout, Util, Templates) {

    var detail, $search, $togglePinned,
        $toggleCluster, $clusterSelection, selectedCluster,
        $toggleTimespan, $timespanSelection,
        loadThumbsXhr;

    function selectCluster(selection) {
        selection = parseInt(selection);

        // Unselect current
        $clusterSelection.removeClass('cluster-' + selectedCluster);

        // Select 'none'
        if (selection == 0) {
            $clusterSelection.hide();
            $toggleCluster.find('.cluster-icon').show();
        // Select color
        } else if (0 < selection && selection <= 4) {
            $clusterSelection
                .css({display: 'inline-block'})
                .addClass('cluster-' + (selectedCluster - 1));
            $toggleCluster.find('.cluster-icon').hide();
        } else {
            throw 'Unexpected cluster selection';
        }
        // Store value
        $clusterSelection.data('value', selection);
        selectedCluster = selection;
    }

    function selectTimespan(selection) {
        var value = (function(s) {
            switch (s) {
                case 7 : return 'W';
                case 30 : return '1m';
                case 60 : return '2m';
            }
            throw 'Unexpected timespan selection';
        })(parseInt(selection));

        $timespanSelection.html(value).data('value', selection);
    }

    function getControlsStates() {
        var filters = [];
        $('.filter-group').each(function() {
            var $group = $(this),
                values = [];
            $group.find('li a.active').each(function() {
                values.push($(this).data('id'));
            });
            if (values.length) {
                filters.push({id: $group.data('id'), pk: values});
            }
        });

        return {
            time_span: parseInt($timespanSelection.data('value')),
            cluster: parseInt($clusterSelection.data('value')),
            pin: ($('#toggle-pinned').hasClass('active')),
            search_query: $('#search-input').val(),
            filters: filters
        };
    }

    function updateFiltersVisibility($group) {
        var $filters = $group.find('ul a'),
            index = $group.data('gt-index') || 10;
        if (index > $filters.length) {
            $filters.show();
            $group.find('.filter-more').hide();
        } else {
            $filters.filter(':lt(' + index + ')').show();
            $filters.filter(':gt(' + index + ')').hide();
            $group.find('.filter-more').show();
        }
    }

    function updateControlsStates(data) {
        // Cluster
        if (data.hasOwnProperty('cluster')) {
            selectCluster(data['cluster'] || 0); // TODO remove 'or zero' (check server side)
        }
        // Time span
        if (data.hasOwnProperty('time_span')) {
            selectTimespan(data['time_span']);
        }
        // Pinned
        if (data.hasOwnProperty('pin')) {
            Util.toggleClass($togglePinned, 'active', data['pin'])
        }
        // Search
        if (data.hasOwnProperty('search_query')) {
            $('#search-input').val(data['search_query']);
        }
        // Filters
        if (data.hasOwnProperty('filters')) {
            $('#filter-flap').html(Templates.filters.render({groups: data['filters']}));
        }

        $('.filter-group').each(function(i, group) {
            var $group = $(group);
            if (i == 0) {
                $group.find('.filter-toggle').addClass('active');
                $group.find('.filter-filters').addClass('in');
            }
            updateFiltersVisibility($group);
        });
    }

    function loadThumbs() {
        if (loadThumbsXhr) {
            loadThumbsXhr.abort();
        }

        loadThumbsXhr = $.ajax({
            url: '/feed/stream2/filter',
            data: {'data': JSON.stringify(getControlsStates())},
            dataType: 'xml',
            method: 'GET'
        });
        loadThumbsXhr.success(function(xml) {
            var data = $(xml).find('data');
            if (data.length) {
                updateControlsStates(JSON.parse(data.text()));
            }

            var list = $(xml).find('thumb-list');
            if (list.length) {
                $('#list .endless_data').html($(list.text()));
            }

            loadThumbsXhr = null;
        });
        loadThumbsXhr.error(function() {
            console.log('Fail to load thumbs.');
        });
    }

    $(function() {

        $search = $('#search');
        $togglePinned = $('#toggle-pinned');

        $toggleCluster = $('#toggle-cluster');
        $clusterSelection = $('#cluster-selection');

        $toggleTimespan = $('#toggle-timespan');
        $timespanSelection = $('#timespan-selection');


        /** ---------------------------------------------------------------------------------------------
         *  SEARCH BAR
         */
        $('#toggle-search, #close-search').on('click', function() {
            Util.toggleClass($search, 'opened');
        });

        var searchKeyUpTimeout;
        $('#search-input')
            .on('focus', function(e) {
                e.stopPropagation();
                $(e.delegateTarget).parents('form').addClass('active');
            })
            .on('blur', function(e) {
                e.stopPropagation();
                $(e.delegateTarget).parents('form').removeClass('active');
            })
            .on('keyup', function(e) {
                if (searchKeyUpTimeout) {
                    clearTimeout(searchKeyUpTimeout);
                }
                var code = e.keyCode || e.which;
                if (code == 13) { // Enter pressed
                    loadThumbs();
                } else {
                    searchKeyUpTimeout = setTimeout(function() {
                        loadThumbs();
                    }, 1000);
                }
            });


        /** ---------------------------------------------------------------------------------------------
         *  SELECTORS
         */
        $('.selector button').on('click', function(e) {
            //e.stopPropagation();
            var $selector = $(e.target).parents('.selector').eq(0);
            Util.toggleClass($selector, 'opened');
        });

        // Cluster
        $('#cluster').on('click', '.choices .cluster-value', function(e) {
            selectCluster($(e.target).closest('.cluster-value').data('cluster'));
            $toggleCluster.trigger('click');
            loadThumbs();
        });

        // Timespan
        $('#timespan').on('click', '.choices a', function(e) {
            selectTimespan($(e.target).closest('a').data('timespan'));
            $toggleTimespan.trigger('click');
            loadThumbs();
        });

        // Close selectors on click out
        $(window).on('click', function(e) {
            if (0 == $(e.target).closest('.selector').length) {
                $('.selector').removeClass('opened');
            }
        });

        // Toggle pinned
        $togglePinned.on('click', function(e) {
            Util.toggleClass($(e.delegateTarget), 'active');
            loadThumbs();
        });


        /** ---------------------------------------------------------------------------------------------
         *  FILTERS
         */
        $('#filter-flap')
            .on('click', '.filter-toggle', function(e) {
                var $group = $(e.target).parents('.filter-group').eq(0);
                if (Util.toggleClass($group, 'active')) {
                    $group.find('.filter-filters').collapse('show');
                } else {
                    $group.data('gt-index', 10)
                        .find('.filter-filters').collapse('hide');
                    updateFiltersVisibility($group);
                }
            })
            .on('click', '.filter-group ul a', function(e) {
                Util.toggleClass($(e.target), 'active');
                loadThumbs();
            })
            .on('click', '.filter-more', function(e) {
                var $group = $(e.target).closest('.filter-group');
                var index = $group.data('gt-index') || 10;
                $group.data('gt-index', index + 10);
                updateFiltersVisibility($group);
            });


        /** ---------------------------------------------------------------------------------------------
         *  THUMBS
         */
        $('#list, #detail')
            .on('click', '.thumb-pin', function(e) {
                e.stopPropagation();

                var $button = $(this);
                var $thumb = $(e.target).parents('.thumb').eq(0);

                var pinXhr = $.ajax({
                    type: 'POST',
                    url: '/user/paper/pin',
                    data: {'pk': $thumb.data('id'), 'source': window.location.pathname}
                });
                pinXhr.success(function(json) {
                    if (json.hasOwnProperty('is_liked')) {
                        $button.toggleClass('active', json['is_liked']);
                    }
                });
                pinXhr.fail(function() {
                    console.log('ERROR');
                });

                return false;
            })
            .on('click', '.thumb-remove', function(e) {
                e.stopPropagation();

                var $thumb = $(e.target).parents('.thumb').eq(0);

                var pinXhr = $.ajax({
                    type: 'POST',
                    url: '/user/paper/ban',
                    data: {'pk': $thumb.data('id'), 'source': window.location.pathname}
                });
                pinXhr.success(function(json) {
                    if (json.hasOwnProperty('is_liked') && json['is_ticked']) {
                        $thumb.remove();
                    }
                });
                pinXhr.fail(function() {
                    console.log('ERROR');
                });

                return false;
            });

        // Endless scroll
        $.endlessPaginate({
            paginateOnScroll: true,
            paginateOnScrollMargin: 10,
            onClick: function (context) {
                context.extraData = JSON.stringify(getControlsStates());
                context.url = window.location.href;
            },
            onCompleted: function (fragment) {
            }
        });


        /** ---------------------------------------------------------------------------------------------
         *  DETAIL
         */
        detail = new Detail();
        $(detail)
            .on('etalia.detail.loading', function() {
                Layout.setBusy();
            })
            .on('etalia.detail.loaded', function() {
                Layout.setAvailable();
            });
        $('.document').on('click', '.thumb .title a', function(e) {
            e.preventDefault();

            detail.load($(e.target).closest('.thumb'));

            return false;
        });

    });
});

