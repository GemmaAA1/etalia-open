define(['jquery', 'app/util/templates', 'app/util/utils'], function ($, Templates, Utils) {

    var statesRegistry = {
        id: [],
        opened: [],
        count: []
    };
    statesRegistry.has = function(id) {
        return 0 <= this.id.indexOf(id);
    };
    statesRegistry.add = function(id, opened, count) {
        if (this.has(id)) {
            this.set(id, opened, count);
        }
        this.id.push(id);
        this.opened.push(opened);
        this.count.push(count);
    };
    statesRegistry.set = function(id, opened, count) {
        if (!this.has(id)) {
            this.add(id, opened, count);
        }
        var index = this.id.indexOf(id);
        this.opened[index] = opened;
        this.count[index] = count;
    };
    statesRegistry.get = function(id) {
        if (!this.has(id)) {
            throw 'Undefined group "' + id + '"';
        }
        var index = this.id.indexOf(id);
        return {
            opened: this.opened[index],
            count: this.count[index]
        }
    };
    statesRegistry.setOpened = function(id, opened) {
        if (!this.has(id)) {
            throw 'Undefined group "' + id + '"';
        }
        var index = this.id.indexOf(id);
        this.opened[index] = opened;
    };
    statesRegistry.getOpened = function(id) {
        if (!this.has(id)) {
            throw 'Undefined group "' + id + '"';
        }
        return this.opened[this.id.indexOf(id)];
    };
    statesRegistry.setCount = function(id, count) {
        if (!this.has(id)) {
            throw 'Undefined group "' + id + '"';
        }
        var index = this.id.indexOf(id);
        this.count[index] = count;
    };
    statesRegistry.getCount = function(id) {
        if (!this.has(id)) {
            throw 'Undefined group "' + id + '"';
        }
        return this.count[this.id.indexOf(id)];
    };


    var Filters = function (options) {
        this.config = $.extend({
            debug: false,
            element: '#filter-flap'
        }, options);

        this.$element = $(this.config.element);
        this.$toggle = $(this.config.selection);
        this.$selection = $(this.config.selection);

        this.selected = null;
    };

    Filters.prototype.init = function() {
        var that = this,
            $body = $('body');

        // Initial groups states
        this.$element.find('.filter-group').each(function(i, group) {
            var $group = $(group),
                count = $group.find('.filter-filters ul a:visible').length;

            statesRegistry.add($group.data('id'), $group.hasClass('active'), count);
            that.applyUserSelection($group);
        });

        this.$element
            .on('click', '.filter-toggle', function(e) {
                var $group = $(e.target).closest('.filter-group'),
                    groupId = $group.data('id');

                if (Utils.toggleClass($group, 'active')) {
                    $group.find('.filter-filters').collapse('show');
                    statesRegistry.setOpened(groupId, true);
                } else {
                    $group.find('.filter-filters').collapse('hide');
                    statesRegistry.setOpened(groupId, false);

                    //updateFiltersVisibility($group);
                }
            })
            .on('click', '.filter-group ul a', function(e) {
                var $a = $(e.target).closest('a'),
                    $group = $a.closest('.filter-group'),
                    active = Utils.toggleClass($a, 'active');

                $body.trigger('etalia.control.filters.change', {
                    value:  $a.data('id'),
                    label:  $a.attr('title'),
                    group:  $group.data('id'),
                    active: active
                });
            })
            .on('click', '.filter-more', function(e) {
                var $group = $(e.target).closest('.filter-group'),
                    groupId = $group.data('id'),
                    count = statesRegistry.getCount(groupId);

                statesRegistry.setCount(groupId, count + 10);

                that.applyUserSelection($group);
            });

        return this;
    };

    Filters.prototype.applyUserSelection = function($group) {
        var $filters = $group.find('ul li'),
            states = statesRegistry.get($group.data('id'));

        if (states.opened) {
            $group.addClass('active')
                .find('.filter-filters').addClass('in');
        } else {
            $group.removeClass('active')
                .find('.filter-filters').removeClass('in');
        }

        if (states.count >= $filters.length) {
            $filters.show();
            $group.find('.filter-more').hide();
        } else {
            $filters.filter(':lt(' + states.count + ')').show();
            $filters.filter(':gt(' + (states.count-1) + ')').hide();
            $group.find('.filter-more').show();
        }

        return this;
    };

    Filters.prototype.render = function(groups) {
        var that = this;

        this.$element
            .html(Templates.filters.render({groups: groups}))
            .find('.filter-group')
                .each(function(i, group) {
                    that.applyUserSelection($(group));
                });

        return this;
    };

    Filters.prototype.setValue = function() {
        throw 'Not implemented';
    };

    Filters.prototype.getValue = function() {
        var value = [];
        this.$element.find('.filter-group').each(function(i, group) {
            var $group = $(group),
                keys = [];
            $group.find('li a.active').each(function() {
                keys.push($(this).data('id'));
            });
            if (keys.length) {
                value.push({id: $group.data('id'), pk: keys});
            }
        });
        return value;
    };

    return Filters;
});
