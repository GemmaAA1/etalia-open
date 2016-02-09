define(
    ['jquery', 'app/api', 'app/ui/detail', 'app/ui/layout', 'app/util/utils', 'app/templates', 'endless', 'bootstrap'],
    function($, Api, Detail, Layout, Util, Templates) {

    var $body, detail, $search, $togglePinned,
        $toggleCluster, $clusterSelection, selectedCluster,
        $toggleTimespan, $timespanSelection,
        openedFiltersGroups,
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
                case 7 :  return 'W';
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
        // List title
        if (data.hasOwnProperty('number_of_papers')) {
            $('.list-title span').html(data['number_of_papers']);
        }
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
            if (0 <= openedFiltersGroups.indexOf($group.data('id'))) {
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
            url: $('#list').data('load-url'),
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
            console.error('Thumbs request failed');
        });
    }

    $(function() {
        $body = $('body');
        $search = $('#search');
        $togglePinned = $('#toggle-pinned');

        $toggleCluster = $('#toggle-cluster');
        $clusterSelection = $('#cluster-selection');

        $toggleTimespan = $('#toggle-timespan');
        $timespanSelection = $('#timespan-selection');

        openedFiltersGroups = [];
        $('.filter-group.active').each(function() {
            openedFiltersGroups.push($(this).data('id'));
        });

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
                var $group = $(e.target).parents('.filter-group').eq(0),
                    groupId = $group.data('id');
                if (Util.toggleClass($group, 'active')) {
                    $group.find('.filter-filters').collapse('show');
                    if (0 > openedFiltersGroups.indexOf(groupId)) {
                        openedFiltersGroups.push(groupId);
                    }
                } else {
                    $group.data('gt-index', 10)
                        .find('.filter-filters').collapse('hide');
                    var index = openedFiltersGroups.indexOf(groupId);
                    if (0 <= index) {
                        openedFiltersGroups.splice(index, 1);
                    }
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
                var $thumb = $(e.target).parents('.thumb').eq(0);

                Api.pin($thumb.data('id'), window.location.pathname);

                e.stopPropagation();
                return false;
            })
            .on('click', '.thumb-ban', function(e) {
                var $thumb = $(e.target).parents('.thumb').eq(0);

                Api.ban($thumb.data('id'), window.location.pathname);

                e.stopPropagation();
                return false;
            });

        // Endless scroll
        $.endlessPaginate({
            paginateOnScroll: true,
            paginateOnScrollMargin: 10,
            onClick: function (context) {
                context.extraData = JSON.stringify(getControlsStates());
                context.urlOrigin = window.location.href;
            }
        });


        /** ---------------------------------------------------------------------------------------------
         *  DETAIL
         */
        detail = new Detail();
        detail.init();

        $('.document').on('click', '.thumb .title a', function(e) {
            detail.load($(e.target).closest('.thumb'));

            e.preventDefault();
            return false;
        });

        // Events
        $body
            .on('etalia.detail.loading', function() {
                Layout.setBusy();
            })
            .on('etalia.detail.loaded', function() {
                Layout.setAvailable();
            });
    });
});

