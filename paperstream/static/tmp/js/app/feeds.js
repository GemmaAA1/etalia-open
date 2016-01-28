define(
    ['jquery', 'app/ui/detail', 'app/util/utils', 'app/ui/layout', 'bootstrap'],
    function($, Detail, Util) {

    var detail, $search,
        $toggleCluster, $clusterSelection, selectedCluster,
        $toggleTimespan, $timespanSelection;

    $(function() {

        detail = new Detail();

        $search = $('#search');

        $toggleCluster = $('#toggle-cluster');
        $clusterSelection = $('#cluster-selection');

        $toggleTimespan = $('#toggle-timespan');
        $timespanSelection = $('#timespan-selection');

        // Toggle pinned
        $('#toggle-pinned').on('click', function(e) {
            Util.toggleClass($(e.delegateTarget), 'active');
        });

        // Toggle search bar
        $('#toggle-search, #close-search').on('click', function() {
            Util.toggleClass($search, 'opened');
        });

        // Search bar active
        $('#search-input')
            .on('focus', function(e) {
                e.stopPropagation();
                $(e.delegateTarget).parents('form').addClass('active');
            })
            .on('blur', function(e) {
                e.stopPropagation();
                $(e.delegateTarget).parents('form').removeClass('active');
            });

        // Toggle selector
        $('.selector button').on('click', function(e) {
            //e.stopPropagation();
            var $selector = $(e.target).parents('.selector').eq(0);
            Util.toggleClass($selector, 'opened');
        });

        // Cluster selection
        $('#cluster').on('click', '.choices .cluster-value', function(e) {
            //e.stopPropagation();
            var $target = $(e.target);
            if (!$target.hasClass('cluster-value')) {
                $target = $target.parents('.cluster-value').eq(0);
            }

            var selection = $target.data('cluster');
            $clusterSelection.removeClass('cluster-' + selectedCluster);
            selectedCluster = selection;
            if (selectedCluster == 'none') {
                $clusterSelection.hide();
                $toggleCluster.find('.cluster-icon').show();
            } else {
                $clusterSelection.css({display: 'inline-block'}).addClass('cluster-' + selectedCluster);
                $toggleCluster.find('.cluster-icon').hide();
            }

            $toggleCluster.trigger('click');
        });

        // Timespan selection
        $('#timespan').on('click', '.choices a', function(e) {
            //e.stopPropagation();
            var value = (function(selection) {
                switch (selection) {
                    case 0 : return 'L';
                    case 1 : return 'W';
                    case 2 : return '1m';
                    case 3 : return '2m';
                }
            })($(e.target).closest('a').data('timespan'));

            $timespanSelection.html(value);
            $toggleTimespan.trigger('click');
        });

        // Close selectors on click out
        $(window).on('click', function(e) {
            if (0 == $(e.target).closest('.selector').length) {
                $('.selector').removeClass('opened');
            }
        });


        // Filters
        $('.filter-toggle').on('click', function(e) {
            var $group = $(e.delegateTarget).parents('.filter-group').eq(0);
            if (Util.toggleClass($group, 'active')) {
                $group.find('.filter-filters').collapse('show');
            } else {
                $group.find('.filter-filters').collapse('hide');
            }
        });
        $('.filter-group ul a').on('click', function(e) {
            Util.toggleClass($(e.delegateTarget), 'active');
        });
        $('.filter-group.active .collapse').collapse('show');


        /* TODO MOVE */
        function getCookie(name) {
            var cookieValue = null;
            if (document.cookie && document.cookie != '') {
                var cookies = document.cookie.split(';');
                for (var i = 0; i < cookies.length; i++) {
                    var cookie = jQuery.trim(cookies[i]);
// Does this cookie string begin with the name we want?
                    if (cookie.substring(0, name.length + 1) == (name + '=')) {
                        cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                        break;
                    }
                }
            }
            return cookieValue;
        }
        var csrftoken = getCookie('csrftoken');
        function csrfSafeMethod(method) {
            // these HTTP methods do not require CSRF protection
            return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
        }
        $.ajaxSetup({
            beforeSend: function (xhr, settings) {
                if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
                    xhr.setRequestHeader("X-CSRFToken", csrftoken);
                }
            }
        });

        // Thumbs
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
                pinXhr.done(function(json) {
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

                var $button = $(this);
                var $thumb = $(e.target).parents('.thumb').eq(0);

                var pinXhr = $.ajax({
                    type: 'POST',
                    url: '/user/paper/ban',
                    data: {'pk': $thumb.data('id'), 'source': window.location.pathname}
                });
                pinXhr.done(function(json) {
                    if (json.hasOwnProperty('is_liked') && json['is_ticked']) {
                        $thumb.remove();
                    }
                });
                pinXhr.fail(function() {
                    console.log('ERROR');
                });

                return false;
            });


        // Detail
        $('.document').on('click', '.thumb .title a', function(e) {
            e.preventDefault();
            detail.load();
            return false;
        });
        $('#detail-close, #backdrop').on('click', function() {
            detail.close();
        });
        $('#detail-next, #detail-prev').on('click', function() {
            detail.load();
        });
    });
});

